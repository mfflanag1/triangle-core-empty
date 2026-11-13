while [ $(wc -l < ip40.jsonl) -lt 20 ]; do sleep 5; done
pkill -f "sim.py 40 6"
/Users/maxf/miniforge3/bin/python3 sim.py 48 6 12 0 300 1 ip48.jsonl
