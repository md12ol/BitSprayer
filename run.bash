#!/bin/bash

: '
Parameters are:
States: [8. 12, 16]
Population: [36, 72]
Mutations: [2, 3, 4]
'

>nohup.out

#nohup ./cmake-build-rl1/DoorTest 16 36 3 2 &
#nohup ./cmake-build-rl1/DoorTest 16 36 3 3 &
#nohup ./cmake-build-rl1/DoorTest 16 36 4 3 &

#nohup ./cmake-build-rl2/DoorTest 16 36 3 2 &
#nohup ./cmake-build-rl2/DoorTest 16 36 3 3 &
#nohup ./cmake-build-rl2/DoorTest 16 36 4 3 &

#nohup ./cmake-build-wsl/DoorTest 12 72 4 2 &
#nohup ./cmake-build-wsl/DoorTest 16 72 2 2 &
#nohup ./cmake-build-wsl/DoorTest 16 72 3 2 &
#nohup ./cmake-build-wsl/DoorTest 16 72 4 2 &
#nohup ./cmake-build-wsl/DoorTest 12 72 3 2 &
#nohup ./cmake-build-wsl/DoorTest 12 72 3 3 &

nohup ./cmake-build-skybytes/DoorTest 16 36 3 2 &
#nohup ./cmake-build-skybytes/DoorTest 16 36 3 3 &

#nohup ./cmake-build-skybytes2/DoorTest 16 36 3 2 &
#nohup ./cmake-build-skybytes2/DoorTest 16 36 3 3 &
