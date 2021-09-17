## Usage

First, make sure you have [`sshpass`](https://gist.github.com/arunoda/7790979) installed.

Next, set the following environment variables, which are used to log into the LBL server.
* `LBL_USER`
* `LBL_PASSWORD`

Optionally, set up a Python virtual environment and install the requirements (or just install the requirements):
```
python3 -m venv _venv
source _venv/bin/activate
python3 -m pip install -r requirements.txt
```

Before running, you will need to mirror some of the files from the LBL server to your local machine:
```
make prepare
```

Now you're ready to go! Running the following will create `build/alignments.csv`:
```
make align
```
