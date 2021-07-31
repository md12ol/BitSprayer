#!/bin/bash

: '
Parameters are:
States: [8. 12, 16]
Population: [36, 72]
Mutations: [2, 3, 4]
'

nohup ./DoorTest 8 36 2 &
nohup ./DoorTest 8 36 3 &
nohup ./DoorTest 8 36 4 &
nohup ./DoorTest 12 36 2 &
nohup ./DoorTest 12 36 3 &
nohup ./DoorTest 12 36 4 &
nohup ./DoorTest 16 36 2 &
nohup ./DoorTest 16 36 3 &
nohup ./DoorTest 16 36 4 &
