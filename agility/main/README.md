# Agility Robotics Cassie
This software stack contains the core real-time and simulation tools for developing on a Cassie robot.
 
## Requirements
* External Repositories
  * agility-matlab-toolkit
* MATLAB 2016b
  * MATLAB 9.1
  * Simulink 8.8
  * Simscape 4.1
  * Simulink Coder 8.11
  * Simulink Real-Time 6.5
  * (Optional) Symbolic Math Toolbox 7.1
* Supported Compiler (https://www.mathworks.com/support/compilers.html)
  * Microsoft Visual C++ 2015 Professional

## Usage
The repository is managed using Simulink Project with git integration. Open MATLAB and navigate to the **agility-cassie** repository, then open the Simulink Project
```matlab
open AgilityCassie.prj
```
This will start Simulink Project and run a startup script that will set all neccesary paths as well as create build directories for code generation. If this is a fresh install on a computer you will also have to add the Target PC settings and open the Simulink Real-Time explorer
```matlab
addTargetPc
slrtexplr
```
Once added, right-click on the **CassieV3** Target PC and click **Set as Default Target**. The unused **TargetPC** can be deleted for cleaniness.

## Testing
Verify the simulation enviroment works running the **NullSimulator.slx** model.
```matlab
sim NullSimulator
```
You should see an animation of the robot falling through the ground in the Mechanics Explorer followed by a summary of the real-time performance.

Next, verify the real-time environment works by building the **NullRealTime.slx** model.
```matlab
slbuild NullRealTime
 ```
 Once complete, there should be three files in the **Build -> Code -> NullRealTime_emb** folder. These are the standalone real-time kernel and application files for booting the robot from a USB.
