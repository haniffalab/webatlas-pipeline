Tested with Python 3.8

From the project root directory:

```
python -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install -r ./docker/requirements.txt 
cd sphinx/
make html
```