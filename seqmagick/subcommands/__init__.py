commands = 'convert', 'info', 'mogrify', 'muscle'

def itermodules(root=__name__):
    for command in commands:
        yield command, __import__('%s.%s' % (root, command), fromlist=[command])
