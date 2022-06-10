import functools
import warnings


def call_on_import(func):
    """Declare a decorator that will run `func` when imported."""
    func()


def get_package_file_path(from_package, relative_path):
    """Use source of a python package to locate and cache the address of a file."""
    from pkg_resources import resource_filename

    return resource_filename(from_package, relative_path)


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
