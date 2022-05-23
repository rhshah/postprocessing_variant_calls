# Main Import 
import typer
import vardict.vardict_process
import logging
import time
# setup logger 
logger = logging.getLogger("filter")
app = typer.Typer()


# Vardict filter 
app.add_typer(vardict.vardict_process.app, name="vardict")
# Haplo filter 
#### Insert Here ####

# Main App 
if __name__ == "__main__":
    app()
