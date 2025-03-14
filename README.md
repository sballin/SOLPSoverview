# SOLPSoverview

## setup

    git clone https://github.com/sballin/SOLPSoverview.git
    cd SOLPSoverview

    # Create a virtual environment in a new folder "venv"
    python -m venv venv
    
    # Activate the virtual environment
    source venv/bin/activate.csh

    # Install required packages
    pip install -r requirements.txt

Add the commands to your path in ~/.tcshrc (assuming you cloned in ~/SOLPSoverview):

    set path = (~/SOLPSoverview/bin $path)

## commands

In your SOLPS run directory, you can do...

### plotall

generate pdf with all plots

### run

    run 10

will run the case as set up 10 times and create directories history/1, ..., history/10 containing the full archive of SOLPS files, pdf plots, .h5 file with data for monitoring webpage

### quicksave

generate h5 file with data for monitoring webpage

### post

do plotall and quicksave