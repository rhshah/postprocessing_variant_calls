# Contributing to Post-Processing 

1. Follow `README.md` for package setup. 
2. Determine the file type that you'd like to add a post-processing operation. You can run `pv --help` to see supported file types. 

3. If the file type exists, for example maf, add the operation as a command for that file type. You can run `pv <filetype> --help` to see available commands.

4. If the command doesn't require sub-commands, for example `pv maf <command>`, you can simply add a new command block to `postprocesssing_variant_calls/<filetype>/main.py`:
    ```
    #TODO change command name and help
    @app.command(
        "germline_status", help="Tag a variant in a MAF file as germline based on VAF value"
    )
    #TODO change command function name 
    #TODO change argument names
    def germline_status(
        maf: Path = typer.Option(
            ...,
            "--maf",
            "-m",
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            help="MAF file to subset",
        ),
        output_maf: Path = typer.Option(
            "output.maf", "--output", "-o", help="Maf output file name."
        ),
        separator: str = typer.Option(
            "tsv",
            "--separator",
            "-sep",
            help="Specify a seperator for delimited data.",
            callback=check_separator,
        ),
    ):
        #TODO change file class to 
        # prep maf
        mafa = MAFFile(maf, separator)
        #TODO change file type operation 
        #TODO you will need to add the operation(s) as a method to the file's class
        mafa = mafa.tag("germline_status")
        #TODO change file writing to handle expected type
        mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
        return 0

    ```
    After adding the block you will need to resolve the `TODOs` in the code block. 

    The most important of these is adding your operation as a method to the file's class. 

    The file's class can be found in: `postprocesssing_variant_calls/<filetype>/helper.py`. 

    If you are un-familiar with classes and methods you may need to read: https://www.w3schools.com/python/python_classes.asp. 

5. If the command does require sub-commands, for example `pv maf <command> <sub-command>`, you must add a new process script under a new sub-command folder. 

    To do this run: 
    ```
    mkdir postprocessing_variant_calls/<filetype>/<command>
    touch postprocessing_variant_calls/<filetype>/<command>/__init__.py
    cp postprocessing_variant_calls/data/method_process.py postprocessing_variant_calls/<filetype>/<command>/<command>_process.py
    ```

    Then, create your sub-commands by completing the `TODOs` noted in `method_process.py``. 

    There also may need to be additional changes depending on the requirements of the new operation such as adding helper functions. These can be added to `postprocessing_variant_calls/<filetype>/helper.py`

    Lastly, you will need to add this block to the top of  `postprocessing_variant_calls/<filetype>/main.py`: 

    ```
    from .<command> import <command>_process
    app.add_typer(<command>_process.app, name="<command>", help="<command_help>")
    ```

    Note: the new command is specified by the block added to the `postprocessing_variant_calls/<filetype>/main.py`. While the sub-commands are specified in the copy of `method_process.py`. 

6. If the file type does not exist, you must create a new class for the file type. Please reach out to @buehlere if this type of feature is needed.
