# Project context database

`project_context.sqlite` holds the Emory Astro Research project's conversation and docs, so any thread can query them instead of re-reading chats.

## Tables
- `threads`: one row per thread (title, status, last activity).
- `messages`: every project-chat and thread message. `body` is raw; `text` is the readable text (status checklists and other structured posts are unwrapped); `kind` is set for structured posts (status, thread_started, skill_proposal, ...). `in_timeline=1` means it appeared in the project chat.
- `documents`: repo `*.md`/`*.txt` docs, keyed `branch:path` (`main:` and the research branch `claude/project-thread-alvlwx:`), skipping data/outputs/results.
- `search`: FTS5 full-text index over messages and docs.
- `sync_state`: per-thread fetch cursors and `last_sync_at`.

## Use
```
python3 ingest.py --stats
python3 ingest.py --search 'J1604 AND W4'
sqlite3 project_context.sqlite "select created_at, author, substr(text,1,120) from messages where thread_id like '%19N33y%' order by created_at desc limit 5"
```

## Sync
An hourly scheduled task wakes the "Project context database sync" thread. It lists threads active since `last_sync_at`, fetches only messages after each stored cursor, writes the fetch JSON to `inbox/` and runs `python3 ingest.py inbox/<file>.json`, re-ingests repo docs, then `--mark-sync`. Ingest is idempotent (upsert by message id / doc hash).

Limits: it is an hourly poll, not a live trigger, so new context lands within about an hour. Edits to messages older than a thread's cursor are not picked up by incremental fetches.
