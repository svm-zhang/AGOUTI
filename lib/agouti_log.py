import logging

class AGOUTI_LOG(object):
	def __init__(self):
		pass

	def create_logger(self, logFile=None):
		pass

class PROGRESS_METER(AGOUTI_LOG):
	def __init__(self, loggerName):
		logLevel = logging.INFO
		self.logger = logging.getLogger(loggerName.upper())
		self.logger.setLevel(logLevel)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		consoleHandler = logging.StreamHandler()
		consoleHandler.setLevel(logLevel)
		consoleHandler.setFormatter(formatter)
		self.logger.addHandler(consoleHandler)

	def add_file_handler(self, logFile, mode='w'):
		fileHandler = logging.FileHandler(logFile, mode=mode)
		fileHandler.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		fileHandler.setFormatter(formatter)
		self.logger.addHandler(fileHandler)

class DEBUG(AGOUTI_LOG):
	def __init__(self, loggerName, logFile, mode='w'):
		logLevel = logging.DEBUG
		logger = logging.getLogger(loggerName.upper())
		logger.setLevel(logLevel)
		fileHandler = logging.FileHandler(logFile, mode=mode)
		fileHandler.setLevel(logLevel)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
		fileHandler.setFormatter(formatter)
		logger.addHandler(fileHandler)
		return logger

	def create_logger(self, logFile):
		logger = logging.getLogger(self.loggerName)
		logger.setLevel(self.logLevel)
		formatter = logging.Formatter('%(name)s\t%(message)s')
		fileHandler = logging.FileHandler(logFile, mode='w')
		fileHandler.setLevel(self.logLevel)
		fileHandler.setFormatter(formatter)
		logger.addHandler(fileHandler)
		return logger
