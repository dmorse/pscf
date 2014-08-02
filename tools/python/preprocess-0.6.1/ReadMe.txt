        ReadMe file for 'preprocess', a multi-language preprocessor.

==== Table of Contents

Introduction
Install Notes
Changelog
Known Bugs
TODO list/Feature Requests


==== Introduction

'preprocess' is a multi-language preprocessor.


==== Install Notes

1. Download the latest 'preprocess' package for your platform from
   http://starship.python.net/~tmick.
2. Unpack the package in a temporary directory.
3. Run 'python setup.py install' in the created preprocess-<version> directory.


==== Changelog

0.6.1:
    - Fix a bug where preprocessor statements were not ignored when not
      emitting. For example the following should _not_ cause an error:
        # #if 0
        # #error womba womba womba
        # #endif
    - Fix a bug where multiple uses of preprocess.preprocess() in the
      same interpreter would erroneously re-use the same list of
      __preprocessedFiles. This could cause false detection of recursive
      #include's.
    - Fix #include, broken in 0.6.0.

0.6.0:
    - substitution: Variables can now replaced with their defined value
      in preprocessed file content. This is turned OFF by default
      because, IMO, substitution should not be done in program strings.
      I need to add lexing for all supported languages before I can do
      *that* properly. Substitution can be turned on with the
      --substitute command-line option or the subst=1 module interface
      option.
    - Add support for preprocessing HTML files.

0.5.0:
    - Add #error, #define, #undef, #ifdef and #ifndef statements.
    - #include statement, -I command line option and 'includePath'
      module interface option to specify an include path
    - Add __FILE__ and __LINE__ default defines.
    - More strict and more helpful error messages:
        - Lines of the form "#else <expr>" and "#endif <expr>" no longer
          match.
        - error messages for illegal #if-block constructs
        - error messages for use of defined(BAR) instead of
          defined('BAR') in expressions
    - New "keep lines" option to output blank lines for skipped content
      lines and preprocessor statement lines (to preserve line numbers
      in the processed file).

0.4.0:
    - Add #elif preprocessor statement.
    - Add defined() built-in, e.g. #if defined('FOO')

0.3.2:
    - Make #if expressions Python code.
    - Change "defines" attribute of preprocess.preprocess().
    - Add -f|--force option to overwrite given output file.

0.2.0:
    - Add content types for C/C++.
    - Better module documentation.
    - You can define *false* vars on the command line now.
    - 'python setup.py install' works.

0.1.0:
    - First release.



==== Known Bugs

- Id: 1
  Title: Nested #if-blocks will not be properly handled.
  Test-Case: test_nested_bug1.py
  Description:
    In the following code:
         #if 0
             #if 1
                 foo
             #endif
         #endif
    the "foo" line _will_ get printed, even though it should not. This is
    because the parsing code just looks at states[-1][0]==EMIT to
    determine if to emit a line. It must instead look at
    states[*][0]==EMIT to determine if can emit a line.
   
- Id: 2
  Title: Substitution (when turned on via -s) will substitute into
      program strings
  Test-Case: test_subst_bug2.py
  Description:
    That is not the ideal behaviour (ideal being defined by (1) what I
    would expect in my program, and (2) what the C preprocessor does).



==== TODO list/Feature Requests

- test 'includePath' optional argument
- test cases for the command line usage
- test cases for -k|--keep-lines
- patchtemplate functionality: look at how the preprocessor.pl below does
  this, perhaps this is an acceptible poorman's version
  (Somewhat have this, see --substitute added in v0.6.0.)
- stealing from http://software.hixie.ch/utilities/unix/preprocessor/
    - Should I really add #elifdef and #elifndef?
    - #include
        - ensure no recursive #include
        - '-I' argument to allow an #include-path
    - Would #filter/#endfilter be useful?
- be more strict:
    - perhaps add -W<arg> option to allow levels of strictness
      - make #undef of undefined vars an error
- perhaps add #pragma to enable turning on of options: e.g.
    # #pragma substitution
    
