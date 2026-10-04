#!/bin/zsh
# Run queued PHOEBE fits, keeping at most MAXJOBS run_phoebe_candidate.py processes alive.
# Queue file: one job per line, "<designation> <tag> [extra args]"; lines are consumed from the top,
# so jobs can be appended while the queue runs. Usage: src/pipeline/phoebe_queue.sh outputs/candidates/queue.txt [MAXJOBS]
cd "$(dirname "$0")/../.."
Q=$1; MAXJOBS=${2:-11}
while true; do
  if [[ ! -s $Q ]]; then sleep 60; [[ -s $Q ]] || break; fi
  while (( $(pgrep -f run_phoebe_candidate.py | wc -l) >= MAXJOBS )); do sleep 30; done
  line=$(head -1 $Q); sed -i '' 1d $Q
  [[ -z $line || $line == \#* ]] && continue
  args=(${=line}); name=$args[1]; tag=$args[2]; extra=(${args[3,-1]})
  [[ $tag == "-" ]] && tag=""
  echo "$(date -u +%FT%TZ) start $name$tag $extra"
  nohup ~/phoebe-env/bin/python src/pipeline/run_phoebe_candidate.py --designation $name ${tag:+--tag $tag} $extra \
    > outputs/candidates/logs/${name}${tag}_stdout.log 2>&1 &
  sleep 5
done
echo "$(date -u +%FT%TZ) queue empty"
