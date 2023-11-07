Tested with Python 3.10

From the project root directory:

```
python -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install -r ./envs/requirements.txt
pip install -r ./envs/dev/requirements.txt
cd sphinx/
make html
```
