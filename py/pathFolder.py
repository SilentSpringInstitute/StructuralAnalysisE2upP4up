from os import listdir, remove, makedirs, path
from shutil import rmtree




def cleanFolder(prin, l_p_save=[]):

    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            p_remove = prin + filin
            if p_remove in l_p_save:
                continue
            else:
                try: remove(p_remove)
                except: rmtree(p_remove)
    return prin




def createFolder(prin, clean=0):

    if not path.exists(prin):
        makedirs(prin)

    if clean == 1:
        cleanFolder(prin)

    return prin