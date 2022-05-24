# Main Import 
import typer
from .vardict import vardict_process
import logging
import time
# setup logger 
logger = logging.getLogger("filter")
app = typer.Typer()

# Vardict filter 
app.add_typer(vardict_process.app, name="vardict")
# Haplo filter 
#### Insert Here ####

# Main App 
if __name__ == "__main__":
    app()
