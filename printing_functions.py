import numpy as np
from datetime import datetime, timezone
from os import getlogin, path

### README
# this file contains functions to support various print commands

from text_options  import text_options

# TODO fix timezone (or note in output that it's UTC)
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
  print(f"{time}", *args, **kwargs)


# TODO look for python method, believe it's already written
def center(input_string):
  '''
  Helper function to center text on a standard laptop screen
  '''
  spacer = "-"
  screen = 76
  center = (screen - (len(input_string)))//2
  return spacer*center + input_string + spacer*center


def attention(input_string):
  '''
  Helper function to print a string in a large font in a central
  position so that it cannot be missed. Blinking text can be added
  if the user enrolls themself by adding their login handle to the if statement. 
  '''
  print(center("THE FINAL STATE MODE IS"))
  if getlogin() == "ballmond":
    center_val = (76 - 3*len(input_string))//2
    s_1 = text_options["bold_italic_blink"] + text_options["green"] + input_string + text_options["reset"]
    s_2 = text_options["bold_italic_blink"] + text_options["yellow"] + input_string + text_options["reset"]
    s_3 = text_options["bold_italic_blink"] + text_options["purple"] + input_string + text_options["reset"]
    s_full = center_val*" " + s_1+s_2+s_3 + center_val*" " 
    print(s_full)
  print(center(input_string))


