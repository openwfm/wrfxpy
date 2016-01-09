
class GribSource:
  """
  Parent class of all grib sources.
  """

  def __init__(self):
    pass

  def vtables(self):
    """
    Returns the vtables that must be set for use with
    this source.
    """
    return {}

