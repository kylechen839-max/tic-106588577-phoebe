#!/usr/bin/env python3
"""Project context database: ingest project chat/thread fetches and repo docs into SQLite.

Usage:
  python3 ingest.py FILE...          # files holding fetch_thread / fetch_project_timeline /
                                     # list_thread_sessions JSON (raw, or wrapped in tool-output text)
  python3 ingest.py --docs REPO_DIR [REF_LABEL]  # ingest *.md / *.txt docs from a checkout
  python3 ingest.py --stats
  python3 ingest.py --search "query"
Idempotent: messages upsert by id, docs upsert by path+sha.
"""
import hashlib, json, os, sqlite3, subprocess, sys, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.environ.get("CONTEXT_DB", os.path.join(HERE, "project_context.sqlite"))

SCHEMA = """
CREATE TABLE IF NOT EXISTS threads(
  thread_id TEXT PRIMARY KEY, title TEXT, session_id TEXT, status_bucket TEXT,
  resolved INTEGER, created_at TEXT, last_activity_at TEXT);
CREATE TABLE IF NOT EXISTS messages(
  id TEXT PRIMARY KEY, thread_id TEXT, created_at TEXT, edited_at TEXT,
  author TEXT, author_id TEXT, author_name TEXT, kind TEXT, body TEXT, text TEXT,
  in_timeline INTEGER DEFAULT 0, ingested_at TEXT);
CREATE INDEX IF NOT EXISTS messages_thread ON messages(thread_id, created_at);
CREATE TABLE IF NOT EXISTS documents(
  path TEXT PRIMARY KEY, source TEXT, git_ref TEXT, sha256 TEXT, content TEXT, ingested_at TEXT);
CREATE TABLE IF NOT EXISTS sync_state(key TEXT PRIMARY KEY, value TEXT);
CREATE VIRTUAL TABLE IF NOT EXISTS search USING fts5(kind, ref, thread_id, created_at, text);
"""

def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

def connect():
    con = sqlite3.connect(DB)
    con.executescript(SCHEMA)
    return con

def extract_json(raw):
    raw = raw.strip()
    if raw.startswith("{"):
        return json.loads(raw)
    for line in raw.splitlines():
        s = line.strip()
        if s.startswith("{"):
            return json.loads(s)
    raise ValueError("no JSON object found")

def readable(body):
    """Structured agent bodies (status, thread_started, ...) are JSON; pull out kind + text."""
    if body and body.lstrip().startswith("{"):
        try:
            d = json.loads(body)
            if isinstance(d, dict):
                return d.get("kind"), d.get("text") or body
        except ValueError:
            pass
    return None, body

def index(con, kind, ref, thread_id, created_at, text):
    con.execute("DELETE FROM search WHERE kind=? AND ref=?", (kind, ref))
    con.execute("INSERT INTO search VALUES(?,?,?,?,?)", (kind, ref, thread_id, created_at, text))

def ingest_messages(con, msgs, timeline):
    n = 0
    for m in msgs:
        kind, text = readable(m.get("body", ""))
        old = con.execute("SELECT body, edited_at FROM messages WHERE id=?", (m["id"],)).fetchone()
        if old and old == (m.get("body"), m.get("edited_at")):
            if timeline:
                con.execute("UPDATE messages SET in_timeline=1 WHERE id=?", (m["id"],))
            continue
        con.execute("""INSERT INTO messages(id,thread_id,created_at,edited_at,author,author_id,author_name,
                       kind,body,text,in_timeline,ingested_at) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)
                       ON CONFLICT(id) DO UPDATE SET body=excluded.body, text=excluded.text,
                       kind=excluded.kind, edited_at=excluded.edited_at,
                       in_timeline=max(in_timeline, excluded.in_timeline), ingested_at=excluded.ingested_at""",
                    (m["id"], m.get("thread_id"), m.get("created_at"), m.get("edited_at"), m.get("author"),
                     m.get("author_id"), m.get("author_name"), kind, m.get("body"), text,
                     1 if timeline else 0, now()))
        index(con, "message", m["id"], m.get("thread_id"), m.get("created_at"), text)
        n += 1
    return n

def ingest_sessions(con, sessions):
    for s in sessions:
        con.execute("""INSERT INTO threads VALUES(?,?,?,?,?,?,?) ON CONFLICT(thread_id) DO UPDATE SET
                       title=excluded.title, session_id=excluded.session_id, status_bucket=excluded.status_bucket,
                       resolved=excluded.resolved, last_activity_at=excluded.last_activity_at""",
                    (s["thread_id"], s.get("title"), s.get("session_id"), s.get("status_bucket"),
                     1 if s.get("resolved") else 0, s.get("created_at"), s.get("last_activity_at")))
    return len(sessions)

def ingest_file(con, path):
    d = extract_json(open(path, encoding="utf-8").read())
    if "sessions" in d:
        return f"{path}: {ingest_sessions(con, d['sessions'])} threads"
    if "messages" in d:
        timeline = "thread_id" not in d  # fetch_thread results carry a top-level thread_id
        n = ingest_messages(con, d["messages"], timeline)
        if d.get("cursor"):
            key = "timeline_cursor" if timeline else f"thread_cursor:{d['thread_id']}"
            con.execute("INSERT OR REPLACE INTO sync_state VALUES(?,?)", (key, d["cursor"]))
        return f"{path}: {n} new/changed messages of {len(d['messages'])}"
    raise ValueError(f"{path}: unrecognised JSON shape")

def ingest_docs(con, repo, ref=None):
    """Docs are keyed 'ref:path' so several branches can coexist."""
    ref = ref or subprocess.run(["git", "-C", repo, "rev-parse", "--abbrev-ref", "HEAD"],
                                capture_output=True, text=True).stdout.strip()
    files = subprocess.run(["git", "-C", repo, "ls-files", "*.md", "*.txt"],
                           capture_output=True, text=True).stdout.split()
    n = 0
    for rel in files:
        if rel.startswith(("data/", "outputs/", "results/")) or rel.startswith("requirements"):
            continue
        content = open(os.path.join(repo, rel), encoding="utf-8", errors="replace").read()
        key = f"{ref}:{rel}"
        sha = hashlib.sha256(content.encode()).hexdigest()
        old = con.execute("SELECT sha256 FROM documents WHERE path=?", (key,)).fetchone()
        if old and old[0] == sha:
            continue
        con.execute("INSERT OR REPLACE INTO documents VALUES(?,?,?,?,?,?)",
                    (key, "repo", ref, sha, content, now()))
        index(con, "document", key, None, None, content)
        n += 1
    return f"docs from {repo}@{ref}: {n} new/changed of {len(files)}"

def stats(con):
    for t in ("threads", "messages", "documents"):
        print(t, con.execute(f"SELECT count(*) FROM {t}").fetchone()[0])
    for row in con.execute("""SELECT coalesce(t.title, m.thread_id), count(*), max(m.created_at)
                              FROM messages m LEFT JOIN threads t USING(thread_id)
                              GROUP BY m.thread_id HAVING count(*)>1 ORDER BY 3"""):
        print("  ", *row)
    for row in con.execute("SELECT * FROM sync_state"):
        print("  state", *row)

def main(argv):
    con = connect()
    if not argv or argv[0] == "--stats":
        stats(con)
    elif argv[0] == "--search":
        for row in con.execute("SELECT kind, ref, created_at, snippet(search, 4, '[', ']', '…', 12) "
                               "FROM search WHERE search MATCH ? ORDER BY rank LIMIT 20", (argv[1],)):
            print(*row, sep=" | ")
    elif argv[0] == "--docs":
        print(ingest_docs(con, argv[1], argv[2] if len(argv) > 2 else None))
    elif argv[0] == "--mark-sync":
        con.execute("INSERT OR REPLACE INTO sync_state VALUES('last_sync_at', ?)", (argv[1] if len(argv) > 1 else now(),))
    else:
        for p in argv:
            print(ingest_file(con, p))
    con.commit()

if __name__ == "__main__":
    main(sys.argv[1:])
