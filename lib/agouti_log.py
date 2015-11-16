import logging

class PROGRESS_METER(object):
	def __init__(self, loggerName):
		logLevel = logging.INFO
		self.logFile = None
		self.logger = logging.getLogger(loggerName.upper()+" PROGRESS")
		self.logger.setLevel(logLevel)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		consoleHandler = logging.StreamHandler()
		consoleHandler.setLevel(logLevel)
		consoleHandler.setFormatter(formatter)
		self.logger.addHandler(consoleHandler)

	def add_file_handler(self, logFile, mode='w'):
		self.logFile = logFile
		fileHandler = logging.FileHandler(self.logFile, mode=mode)
		fileHandler.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		fileHandler.setFormatter(formatter)
		self.logger.addHandler(fileHandler)

class DEBUG(object):
	def __init__(self, debuggerName, logFile, mode='w'):
		logLevel = logging.DEBUG
		self.debugger = logging.getLogger(debuggerName.upper()+" DEBUG")
		self.debugger.setLevel(logLevel)
		fileHandler = logging.FileHandler(logFile, mode=mode)
		fileHandler.setLevel(logLevel)
		formatter = logging.Formatter('%(name)s - %(message)s')
		fileHandler.setFormatter(formatter)
		self.debugger.addHandler(fileHandler)
