""" scripts for moving a specific set of subdirectory from a set of directories to a new destination."""

import os as os
import glob as glob
import shutil as shutil

ROOT_PATH = '/Volumes/Data/HertzLab/'

source_dir = ROOT_PATH + 'ArrayData/'
dest_dir = ROOT_PATH + 'ArrayImages/'

shutil.copytree(source_dir, dest_dir,
                ignore=shutil.ignore_patterns('*.gps', '*.gpr', '*.txt', '*.xls', '*.xlsx', '*.doc', '*.docx', '*.csv', '*.mat', '*.png', '*.eps', '*.pse', '*.db'))

#  remove image directories from source
# img_dirs = glob.glob(source_dir + '/*/*/IMG') + glob.glob(source_dir + '/*/*/*/IMG') + glob.glob(source_dir + '/*/*/*/*/IMG')
# for d in img_dirs:
#     shutil.rmtree(d)


# remove GPR, GPS and FIG directories from dest.
gpr_dirs = glob.glob(dest_dir + '/*/*/GPR') + glob.glob(dest_dir + '/*/*/*/GPR') + glob.glob(dest_dir + '/*/*/*/*/GPR')
gps_dirs = glob.glob(dest_dir + '/*/*/GPS') + glob.glob(dest_dir + '/*/*/*/GPS') + glob.glob(dest_dir + '/*/*/*/*/GPS')
fig_dirs = glob.glob(dest_dir + '/*/*/Figs') + glob.glob(dest_dir + '/*/*/*/Figs') + glob.glob(dest_dir + '/*/*/*/*/Figs')
doc_dirs = glob.glob(dest_dir + '/*/*/docs') + glob.glob(dest_dir + '/*/*/*/docs') + glob.glob(dest_dir + '/*/*/*/*/docs')
fig_dirs1 = glob.glob(dest_dir + '/*/*/figs') + glob.glob(dest_dir + '/*/*/*/figs') + glob.glob(dest_dir + '/*/*/*/*/figs')


all_dirs = gpr_dirs  + gps_dirs + fig_dirs + doc_dirs + fig_dirs1

for d in all_dirs:
    try:
        shutil.rmtree(d)
    except shutil.Error as e:
        print('Directory not removed. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not removed. Error: %s' % e)