import os
import sys
import logging

class AGOUTI_LOG(object):
	def __init__(self, loggerName):
		self.loggerName = loggerName
		self.logLevel = logging.INFO

	def create_logger(self, logFile=None):
		logger = logging.getLogger(self.loggerName)
		logger.setLevel(self.logLevel)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		if logFile is None:
			fileHandler = logging.NullHandler()
		else:
			fileHandler = logging.FileHandler(logFile, mode='w')
			fileHandler.setLevel(self.logLevel)
			fileHandler.setFormatter(formatter)
		consoleHandler = logging.StreamHandler()
		consoleHandler.setLevel(self.logLevel)
		consoleHandler.setFormatter(formatter)
		logger.addHandler(fileHandler)
		logger.addHandler(consoleHandler)
		return logger

class AGOUTI_DEBUG_LOG(AGOUTI_LOG):
	def __init__(self, loggerName):
		super(AGOUTI_DEBUG_LOG, self).__init__(loggerName)
		self.logLevel = logging.DEBUG

	def create_logger(self, logFile):
		logger = logging.getLogger(self.loggerName)
		logger.setLevel(self.logLevel)
		formatter = logging.Formatter('%(name)s\t%(message)s')
		fileHandler = logging.FileHandler(logFile, mode='w')
		fileHandler.setLevel(self.logLevel)
		fileHandler.setFormatter(formatter)
		logger.addHandler(fileHandler)
		return logger
