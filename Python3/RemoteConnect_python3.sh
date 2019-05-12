#!/bin/bash
echo "Remote workers ready"
#ssh euclid12 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
#ssh euclid13 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
#ssh euclid14 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
#ssh euclid15 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
#ssh euclid16 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
ssh euclid17-gpu "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
ssh euclid18-gpu "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
ssh euclid19-gpu "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
ssh euclid20-gpu "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
#ssh euclid11 "cd Masters_project/ProjectDisperse/Python3 && bash SpawnWorkers_python3.sh 10 0" &
