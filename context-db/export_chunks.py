#!/usr/bin/env python3
"""Export project_context.sqlite as RAG-ready chunks (JSONL: id, text, metadata)."""
import json, os, sqlite3, sys
HERE = os.path.dirname(os.path.abspath(__file__))
db = sqlite3.connect(os.path.join(HERE, "project_context.sqlite"))
out = open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "chunks.jsonl"), "w")
MAX = 2000

def split(text):
    parts, cur = [], ""
    for para in text.split("\n\n"):
        if cur and len(cur) + len(para) > MAX:
            parts.append(cur); cur = ""
        while len(para) > MAX:
            parts.append(para[:MAX]); para = para[MAX:]
        cur = f"{cur}\n\n{para}" if cur else para
    if cur.strip():
        parts.append(cur)
    return parts

titles = dict(db.execute("SELECT thread_id, title FROM threads"))
n = 0
for mid, tid, ts, author, name, kind, text in db.execute(
        "SELECT id, thread_id, created_at, author, author_name, kind, text FROM messages ORDER BY created_at"):
    if kind in ("thread_started", "message_routed", "error", "onboarding") or not (text or "").strip():
        continue
    for i, part in enumerate(split(text)):
        meta = {"source": "project_message", "message_id": mid, "thread_id": tid,
                "thread_title": titles.get(tid, "project chat"), "created_at": ts,
                "author": name or author, "kind": kind or "message", "chunk": i}
        out.write(json.dumps({"id": f"{mid}#{i}", "text": part, "metadata": meta}) + "\n"); n += 1
for path, ref, content in db.execute("SELECT path, git_ref, content FROM documents"):
    for i, part in enumerate(split(content)):
        meta = {"source": "repo_doc", "path": path.split(":", 1)[1], "branch": ref, "chunk": i}
        out.write(json.dumps({"id": f"doc:{path}#{i}", "text": part, "metadata": meta}) + "\n"); n += 1
print(n, "chunks")
