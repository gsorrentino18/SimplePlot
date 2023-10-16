### README ###
# text_options is a short dictionary for adding effects to text in terminal output
# "reset" should always be used after another text option so that output after the
# completion of the plotting program remains normal. For example:
# my_string = "i am green and underlined"
# print(text_options["green"] + text_options["uline"] + my_string + text_options["reset"])
# is a valid usage.
# INFO: \033[ is an ANSI escape sequence that usually works for Mac and Linux
# the escape sequence is followed by some number denoting an option, and terminated
# with an m. To string multiple options, see the entry for "bold_italic_blink".

text_options = {
  "reset" : "\033[m",
  "bold"  : "\033[1m",
  "uline" : "\033[4m",
  "blink" : "\033[5m",
  "bold_italic_blink" : "\033[1m\033[4m\033[5m",
  "green"  : "\033[92m",
  "yellow" : "\033[93m",
  "purple" : "\033[94m",
}
