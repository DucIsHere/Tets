import sys
import re
import builtins
from time import sleep

try:
    from rich import print as rprint
except Exception:
    tag_re = re.compile(r'\[/?[^\]]+\]')

    def _strip_rich_tags(text: str) -> str:
        return tag_re.sub('', text)

    def rprint(*args, **kwargs):
        sep = kwargs.pop('sep', ' ')
        end = kwargs.pop('end', '\n')
        file = kwargs.pop('file', None)
        flush = kwargs.pop('flush', False)
        text = sep.join(str(a) for a in args)
        text = _strip_rich_tags(text)
        builtins.print(text, end=end, file=file, flush=flush)

def hide_cursor():
    sys.stdout.write("\033[?25l")
    sys.stdout.flush()

def show_cursor():
    sys.stdout.write("\033[?25h")
    sys.stdout.flush()
    
def printLyrics():
    Lyrics = [
        ("Yo te espero sin na' ", 0.05),
        ("Dime qué es lo que pasa, yo quiero", 0.09),
        ("Pienso quedarme contigo papi", 0.06),
        ("Puedo dar el control contigo papi", 0.065),
        ("Dime dime yo te digo papi,", 0.055),
        ("Ay papi, ay papi", 0.1),
        ("De la fuma mala fuma", 0.05),
        ("Se la pasan diciendo que soy mala", 0.05),
        ("Porque no me aguanto drama de nada", 0.05),
        ("Dime si viense, a mi cama, a mi cama", 0.05),
        
    ]
    
    delays = [0.05, 0.9, 0.08, 0.05, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05]

    hide_cursor()

    try:
        for i, (line, char_delay) in enumerate(Lyrics):
            for char in line:
                if line == "Ay papi, ay papi":
                    rprint(f"[orange4]{char}[/orange4]", end='')
                else:
                    rprint(f"[bold][gold3]{char}[/gold3][/bold]", end='')
                sys.stdout.flush()
                sleep(char_delay)
            rprint()
            sleep(delays[i])
    finally:
        show_cursor()
    
printLyrics()
