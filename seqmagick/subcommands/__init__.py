commands = 'convert', 'info', 'mogrify', 'primer_trim', 'quality_filter', \
        'extract_ids', 'backtrans_align', 'protparam'

def itermodules(root=__name__):
    for command in commands:
        yield (command.replace('_', '-'),
               __import__('%s.%s' % (root, command), fromlist=[command]))
