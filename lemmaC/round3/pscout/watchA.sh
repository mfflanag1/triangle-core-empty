while [ $(grep -c '"q": 256, "K": 8.0' lp.jsonl) -lt 10 ]; do sleep 5; done
pkill -f "sim.py 256 8"
