import functools
import warnings


def call_on_import(func):
    """Declare a decorator that will run `func` when imported."""
    func()


def get_package_file_path(from_package, relative_path):
    """Use source of a python package to locate and cache the address of a file."""
    import importlib_resources

    return str(importlib_resources.files(from_package) / relative_path)


def deprecate_kwargs(deprecated_kwargs=None):
    if deprecated_kwargs is None:
        deprecated_kwargs = set()

    def decorate_deprecate_kwargs(func):
        @functools.wraps(func)
        def wrapper(self_or_cls, *args, **kwargs):
            _deprecate_kwargs(kwargs, deprecated_kwargs)
            return func(self_or_cls, *args, **kwargs)

        return wrapper

    return decorate_deprecate_kwargs


def _deprecate_kwargs(kwargs, deprecated_kwargs):
    added_args = []
    for kwarg in kwargs:
        if kwarg in deprecated_kwargs:
            added_args.append(kwarg)
    if len(added_args) > 1:
        message = (
            "Keyword arguments `{dep_args}` are deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    else:
        message = (
            "Keyword argument `{dep_args}` is deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    if added_args:
        warnings.warn(
            message.format(dep_args=", ".join(added_args)),
            DeprecationWarning,
            3,
        )


def pad_with_wildcards(input_dictionary, max_len, wildcard="*"):
    """Pad empty type or classes with wildcards"""
    types = [f"type{j+1}" for j in range(max_len)]
    classes = [f"class{j+1}" for j in range(max_len)]

    if types[0] in input_dictionary:
        for type_ in types:
            value = input_dictionary[type_]
            if isinstance(value, str) and value.strip() == "":
                input_dictionary[type_] = wildcard

    elif classes[0] in input_dictionary:
        for class_ in classes:
            value = input_dictionary[class_]
            if isinstance(value, str) and value.strip() == "":
                input_dictionary[class_] = wildcard

    return input_dictionary


def get_virtual_ntype_or_nclass_attribs(inputDict):
    """Return classes or types for a virtual type, handling variable length elements."""
    outDict = {"member_classes": [], "member_types": []}
    classes_n = []  # list to sort order by n at end of classn
    types_n = []  # list to sort order by n at end of typen
    for key, val in inputDict.items():
        if key.startswith("class"):
            outDict["member_classes"].append(val)
            classes_n.append(key[5:])  # "class" has 5 letters
        elif key.startswith("type"):
            outDict["member_types"].append(val)
            types_n.append(key[4:])  # "type" has 4 letters
        else:
            outDict[key] = val
    outDict["member_classes"] = [
        x for _, x in sorted(zip(classes_n, outDict["member_classes"]))
    ]  # sort classes
    outDict["member_types"] = [
        x for _, x in sorted(zip(types_n, outDict["member_types"]))
    ]  # sort types
    return outDict
