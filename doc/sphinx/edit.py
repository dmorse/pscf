#!/usr/bin/env python

from FileEditor import *

editor = FileEditor()
editor.setIsTest(False)

#editor.setBlockSize(5)
#editor.setFilter(r'^ *namespace MolMcD *\n *{')
#editor.setOld(r'MolMcD *\n *{( *\n)*')
#editor.setNew(r'Util\n{\n\n')
#for directory in mcmd_directories:
#   editor.editFiles(directory, "*.h")
#   editor.editFiles(directory, "*.cpp")

#editor.setFilter(r'\<p\>')
#editor.setOld(r'\<p\>')
#editor.setNew(r'')
#editor.editFiles(".", "*.rst")
#
#editor.setFilter(r'\<\/p\>')
#editor.setOld(r'\<\/p\>')
#editor.setNew(r'')
#editor.editFiles(".", "*.rst")

editor.setFilter(r'\<li\>')
editor.setOld(r'\<li\>')
editor.setNew(r'\* ')
editor.editFiles(".", "*.rst")

editor.setFilter(r'\<\/li\>')
editor.setOld(r'\<\/li\>')
editor.setNew(r'')
editor.editFiles(".", "*.rst")

