from collections import defaultdict
import functools

# Global registry to store methods
global_methods_registry = defaultdict(list)
global_instances_registry = defaultdict(list)


def register_method(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    print(func, type(func))

    if isinstance(func, (staticmethod, classmethod)):
        class_name = func.__func__.__qualname__.split(".")[0]
    elif "." in func.__qualname__:
        class_name = func.__qualname__.split(".")[0]
    else:
        class_name = None

    # Register the method in the global registry
    global_methods_registry[class_name].append(func)
    if isinstance(func, staticmethod):
        return staticmethod(wrapper)
    elif isinstance(func, classmethod):
        return classmethod(wrapper)
    else:
        return wrapper


def register_class(cls):
    original_init = cls.__init__

    @functools.wraps(original_init)
    def new_init(self, *args, **kwargs):
        class_name = self.__class__.__name__
        global_instances_registry[class_name].append(self)
        original_init(self, *args, **kwargs)

    cls.__init__ = new_init

    # for attr_name, attr_value in cls.__dict__.items():
    #     if callable(attr_value) and not isinstance(attr_value, (staticmethod, classmethod)):
    #         setattr(cls, attr_name, register_method(attr_value))

    return cls


class MethodRunner:

    @staticmethod
    def run_registered_methods():
        for class_name, methods in global_methods_registry.items():
            if class_name is None:
                for method in methods:
                    method()
                continue

            for method in methods:
                if isinstance(method, staticmethod):
                    # Static method, call it directly
                    method.__func__()
                elif isinstance(method, classmethod):
                    class_obj = globals().get(class_name)
                    # Class method, call it with the class
                    method.__func__(class_obj)
                else:
                    for instance in global_instances_registry[class_name]:
                        method(instance)


# Example usage in different classes


@register_method
def basic_method():
    print("Basic Method")


@register_class
class MyClass:

    def __init__(self, val):
        self.val = val

    @register_method
    def method_one(self):
        print(f"MyClass Method One {self.val}")

    @register_method
    @classmethod
    def static_method(cls):
        print("MyClass Static Method")


@register_class
class AnotherClass:
    @register_method
    def another_method(self):
        print("AnotherClass Method")

    @staticmethod
    @register_method
    def another_static_method():
        print("AnotherClass Static Method")


# Testing the functionality
print("Running all registered methods:")
# inst = MyClass(1)
# for method_name in dir(inst):
#     if "__" in method_name:
#         continue
#     method = getattr(inst, method_name, None)
#     if method is not None and callable(method):
#         print(method_name, method, method.__qualname__)
#         print(type(method))
# MethodRunner.run_registered_methods()

inst1 = MyClass(1)
inst2 = MyClass(1)
print(inst2.static_method.__name__)
