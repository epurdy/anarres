"""Miscellaneous useful functions.
"""


NEXT_ID = 0
def get_next_id(typ):
    global NEXT_ID
    rv = '%s_%06d' % (typ, NEXT_ID)
    NEXT_ID += 1
    #print('ID', rv)
    return rv
