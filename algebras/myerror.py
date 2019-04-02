class MyBaseError(Exception):  # TODO: Do not use these exceptions
    pass


class MyClassError(MyBaseError):
    """ This is called when something is wrong in the class hierarchy """
    pass


class MyValueError(MyBaseError):
    pass


class MyTypeError(MyBaseError):
    pass


class MyKeyError(MyBaseError):
    pass
