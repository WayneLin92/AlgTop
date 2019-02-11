class Meta(type):
    pass


class MyClass(metaclass=Meta):
    pass


if __name__ == "__main__":
    a = MyClass()
    print(a)
