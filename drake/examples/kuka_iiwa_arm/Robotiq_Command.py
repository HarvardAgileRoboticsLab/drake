import os
import sys
import lcm

sys.path.append(os.path.abspath("../robotiq-driver"))
import CModelSimpleController
from robo_lcm import lcmt_robotiq_output as outputMsg

class Robotiq_Command(object):
    '''
    Displays simulation time and frequency in the status bar
    using information from the IIWA_STATUS lcm channel.
    '''

    def __init__(self):
        self.command = outputMsg()
        self.lc = lcm.LCM()
        self.command_topic = outputMsg()
        self.command = outputMsg()

    def gripperConnect(self):
      #CModelRtuNode.mainLoop("/dev/ttyUSB0")
      #os.system("CModelRtuNode.py /dev/ttyUSB0")
      CModelSimpleController.publisher()
      pass
    def gripperOpen(self):
      #sendGripperCommand(0, speed=255, force=150)
      self.command = CModelSimpleController.genCommand('o',self.command)
      self.lc.publish("command_topic",self.command.encode())
      #pass
    def gripperClose(self):
      #sendGripperCommand(255, speed=255, force=150)
      self.command = CModelSimpleController.genCommand('c',self.command)
      self.lc.publish("command_topic",self.command.encode())
      #pass
    def gripperReset(self):
      #sendGripperCommand(0, speed=255, force=150)
      #sendGripperCommand('r')
      self.command = CModelSimpleController.genCommand('r',self.command)
      self.lc.publish("command_topic",self.command.encode())
      #pass
    def gripperActivate(self):
      #sendGripperCommand(0, speed=255, force=150)
      #sendGripperCommand('a')
      self.command = CModelSimpleController.genCommand('a',self.command)
      self.lc.publish("command_topic",self.command.encode())