'''
Usage: This program should be launched using the command line specified in the
       kuka_sim.pmd file.
'''

from director import mainwindowapp
from director import robotsystem
from director import applogic
from director import lcmUtils
from PythonQt import QtGui, QtCore
import subprocess
import sys
import os
import lcm
import drake as lcmdrake
#sys.path.append(os.path.abspath("../robotiq-driver"))
print sys.path
#import CModelRtuNode
#import CModelSimpleController
#from robo_lcm import lcmt_robotiq_output as outputMsg
import Robotiq_Command

class KukaSimInfoLabel(object):
    '''
    Displays simulation time and frequency in the status bar
    using information from the IIWA_STATUS lcm channel.
    '''

    def __init__(self, statusBar):
        self.label = QtGui.QLabel('')
        statusBar.addPermanentWidget(self.label)

        self.sub = lcmUtils.addSubscriber('IIWA_STATUS',
                                  lcmdrake.lcmt_iiwa_status, self.onIiwaStatus)
        self.sub.setSpeedLimit(30)

        self.label.text = '[waiting for sim status]'

    def onIiwaStatus(self, msg):
        simTime = msg.utime*1e-6
        simFreq = self.sub.getMessageRate()
        self.label.text = 'Sim freq: %d hz  |  Sim time: %.2f' % (simFreq,
                                                                  simTime)


def makeRobotSystem(view):
    factory = robotsystem.RobotSystemFactory()
    options = factory.getDisabledOptions()
    factory.setDependentOptions(options, usePlannerPublisher=True,
                                         useTeleop=True)
    return factory.construct(view=view, options=options)

def sendGripperCommand(action):

  # create gripper lcm message
  #msg = lcmdrake.lcmt_robotiq_command()
  #msg.utime = int(time.time()*1e6)
  #msg.rACT = 1 # Activate gripper: Gripper only works when this is 0x1 and will complete an auto-callibration after

  #Activates GoTo for the specified position
  #msg.rGTO = 1

  #Automatic release function:Can overide all other commands
  #msg.rATR = 0;

  #Sets the target position for gripper fingers:0x00(0) is open,0xFF(255) is closed
  #msg.rPR = targetPosition;

  #Gripper closing or opening speed:0x00(0) is min speed,0xFF(255) is max speed
  #msg.rSP = speed

  #Gripper final gripping force:0x00(0) is min force,0xFF(255) is max force
  #msg.rFR = force

  #msg.force = force
  #msg.target_position_mm = targetPositionMM

  # send message
  #lcmUtils.publish('SCHUNK_WSG_COMMAND', msg)
  #lcmUtils.publish('ROBOTIQ_COMMAND', msg)
  lc = lcm.LCM()
  command_topic = outputMsg()
  command = outputMsg()
  command = CModelSimpleController.genCommand(action,command) 
  lc.publish("command_topic",command.encode())

def gripperConnect():
  #CModelRtuNode.mainLoop("/dev/ttyUSB0")
  #os.system("CModelRtuNode.py /dev/ttyUSB0")
  #CModelSimpleController.publisher()
  pass
def gripperOpen():
  #sendGripperCommand(0, speed=255, force=150)
  #command = CModelSimpleController.genCommand('o',command)
  #lc.publish("command_topic",command.encode())
  tester.gripperOpen()
  #pass
def gripperClose():
  #sendGripperCommand(255, speed=255, force=150)
  #command = CModelSimpleController.genCommand('c',command)
  #lc.publish("command_topic",command.encode())
  tester.gripperClose()
  #pass
def gripperReset():
  #sendGripperCommand(0, speed=255, force=150)
  #sendGripperCommand('r')
  #command = CModelSimpleController.genCommand('r',command)
  #lc.publish("command_topic",command.encode())
  tester.gripperReset()
  #pass
def gripperActivate():
  #sendGripperCommand(0, speed=255, force=150)
  #sendGripperCommand('a')
  #command = CModelSimpleController.genCommand('a',command)
  #lc.publish("command_topic",command.encode())
  tester.gripperActivate()
  #pass
  
def setupGripperButtons():
  toolBar = applogic.findToolBar('Main Toolbar')
  app.app.addToolBarAction(toolBar, 'Gripper Open', icon='', callback=gripperOpen)
  app.app.addToolBarAction(toolBar, 'Gripper Close', icon='', callback=gripperClose)
  app.app.addToolBarAction(toolBar, 'Gripper Reset', icon='', callback=gripperReset)
  app.app.addToolBarAction(toolBar, 'Gripper Connect', icon='', callback=gripperConnect)
  app.app.addToolBarAction(toolBar, 'Gripper Activate', icon='', callback=gripperActivate)
# create a default mainwindow app
app = mainwindowapp.MainWindowAppFactory().construct()
mainwindowapp.MainWindowPanelFactory().construct(app=app.app, view=app.view)

# load a minimal robot system with ik planning
robotSystem = makeRobotSystem(app.view)

# add the teleop and playback panels to the mainwindow
app.app.addWidgetToDock(robotSystem.teleopPanel.widget,
                        QtCore.Qt.RightDockWidgetArea)
app.app.addWidgetToDock(robotSystem.playbackPanel.widget,
                        QtCore.Qt.BottomDockWidgetArea)

setupGripperButtons()

# show sim time in the status bar
infoLabel = KukaSimInfoLabel(app.mainWindow.statusBar())

# use pydrake ik backend
ikPlanner = robotSystem.ikPlanner
ikPlanner.planningMode = 'pydrake'
ikPlanner.plannerPub._setupLocalServer()

# change the default animation mode of the playback panel
robotSystem.playbackPanel.animateOnExecute = True

# disable pointwise ik by default
robotSystem.ikPlanner.getIkOptions().setProperty('Use pointwise', False)

# set the default camera view
#applogic.resetCamera(viewDirection=[-1,0,0], view=app.view)

tester = Robotiq_Command.Robotiq_Command()

# start!
app.app.start()
