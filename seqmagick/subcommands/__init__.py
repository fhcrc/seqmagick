commands = 'convert', 'info', 'mogrify', 'muscle', 'primer_trim'

def itermodules(root=__name__):
    for command in commands:
        yield (command.replace('_', '-'),
               __import__('%s.%s' % (root, command), fromlist=[command]))
