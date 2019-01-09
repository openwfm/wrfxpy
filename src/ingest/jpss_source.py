#
# Angel Farguell, CU Denver
#

class JPSSError(Exception):
	"""
    Raised when a JPSSSource cannot retrieve JPSS data.
    """
    pass

class JPSSSource(object):
	