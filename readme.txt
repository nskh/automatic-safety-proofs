This walks through steps required to run the examples for our TACAS submission "Automating geometric proofs of collision avoidance via active corners"

On the VM, download and unzip a zip file of the directory:
$ wget https://github.com/nskh/automatic-safety-proofs/archive/clean_artifact.zip
$ unzip clean_artifact.zip

CD into the directory
$ cd automatic-safety-proofs-clean_artifact

Install dependencies
$ chmod +x install_pip.sh
$ ./install_pip.sh

Run examples
$ python3 example_1.py
$ python3 example_2.py
$ python3 example_3.py