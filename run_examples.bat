@echo off
set B=%1
if "%B%"=="" set B=build\Release
if not exist examples\out mkdir examples\out
%B%\flexrotor summary examples\runs.csv --out examples\out
%B%\flexrotor identify examples\runs.csv --setting N2 --out examples\out
%B%\flexrotor landscape examples\runs.csv --setting N2 --params examples\out\params.txt --out examples\out
%B%\flexrotor cv examples\runs.csv --setting N2
%B%\flexrotor orders examples\runs.csv --setting N2 --params examples\out\params.txt
%B%\flexrotor kinematics examples\runs.csv --setting N2 --params examples\out\params.txt
%B%\flexrotor map examples\runs.csv --setting N2 --params examples\out\params.txt --out examples\out
%B%\flexrotor loadmodel examples\runs.csv --out examples\out
%B%\flexrotor pair examples\runs.csv --setting N2 --spacing 0.190 --out examples\out
%B%\flexrotor spring examples\runs.csv --config examples\rig_sggw.cfg --setting N2 --params examples\out\params.txt --out examples\out
%B%\flexrotor virtual examples\runs.csv --config examples\rig_sggw.cfg --setting N2 --params examples\out\params.txt --out examples\out
%B%\flexrotor row examples\runs.csv --config examples\rig_sggw.cfg --setting N2 --params examples\out\params.txt --out examples\out
%B%\flexrotor virtual examples\runs.csv --config examples\rig_sggw.cfg --setting N2 --params examples\out\spring_params_1.0.txt --tag design_1.0 --out examples\out
%B%\flexrotor virtual examples\runs.csv --config examples\rig_sggw.cfg --setting N2 --params examples\out\spring_params_1.5.txt --soil-scale 1.5 --tag design_1.5 --out examples\out
