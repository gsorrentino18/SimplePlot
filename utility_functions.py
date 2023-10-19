import numpy as np
from datetime import datetime, timezone
from os import getlogin, path, makedirs

### README
# this file contains functions to support various print commands
#
# text_options is a short dictionary for adding effects to text in terminal output
# "reset" should always be used after another text option so that output after the
# completion of the plotting program remains normal. For example:
# my_string = "i am green and underlined"
# print(text_options["green"] + text_options["uline"] + my_string + text_options["reset"])
# is a valid usage.
# INFO: \033[ is an ANSI escape sequence that usually works for Mac and Linux
# the escape sequence is followed by some number denoting an option, and terminated
# with an m. To string multiple options, see the entry for "bold_italic_blink".
# more info here: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797

text_options = {
  "reset" : "\033[m",
  "bold"  : "\033[1m",
  "uline" : "\033[4m",
  "blink" : "\033[5m",
  "bold_italic_blink" : "\033[1m\033[4m\033[5m",
  "red"    : "\033[91m",
  "green"  : "\033[92m",
  "yellow" : "\033[93m",
  "purple"   : "\033[94m",
  "pink"   : "\033[95m",
}


def time_print(*args, **kwargs):
  '''
  Helper function to append a time to print statements, the idea being
  to give the user an idea of progress, bottlenecks, and total time of a computation.
  Randomly selected emojis can be added if the user enrolls themself by adding 
  their login handle to the if statement. 
  '''
  emoji = np.random.choice([";) ", ":^)", "<.<", ">.>", ":O ", "^.^", "UwU", "owO"])
  time  = datetime.now(timezone.utc).strftime('%H:%M:%S')
  if getlogin() == "ballmond":
    time = emoji + "  " + time
  print(f"{time} UTC", *args, **kwargs)


def attention(input_string):
  '''
  Helper function to print a string in a large font in a central
  position so that it cannot be missed. Blinking text can be added
  if the user enrolls themself by adding their login handle to the if statement. 
  '''
  screen_width, spacer = 76, "-"
  print("the final state mode is".upper().center(screen_width, spacer))
  if getlogin() == "ballmond":
    # can't use normal center function because escape characters contribute to length of string
    center_val = (screen_width - 3*len(input_string))//2
    s_1 = text_options["bold_italic_blink"] + text_options["green"]  + input_string + text_options["reset"]
    s_2 = text_options["bold_italic_blink"] + text_options["yellow"] + input_string + text_options["reset"]
    s_3 = text_options["bold_italic_blink"] + text_options["purple"] + input_string + text_options["reset"]
    s_full = center_val*" " + s_1+s_2+s_3 + center_val*" " 
    print(s_full)
  print(input_string.center(screen_width, spacer))


def make_directory(directory_name, final_state, testing=False):
  date_and_time  = datetime.now(timezone.utc).strftime('from_%d-%m_at_%H%M')
  directory_name = final_state + "_" + directory_name + "_" + date_and_time
  if testing: directory_name += "_testing"
  if not path.isdir(directory_name):
    makedirs(directory_name)
  else:
    directory_name = "alternate_" + directory_name
    makedirs(directory_name)
    print("WARNING: directory already exists, putting images in alternate: {directory_name}")
  return directory_name
 
