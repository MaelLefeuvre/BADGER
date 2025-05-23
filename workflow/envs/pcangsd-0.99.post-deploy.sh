#!/usr/bin/env bash

set -euxo pipefail

#eval "$(conda shell.bash hook)"
#conda activate "$(basename ${0%.post-deploy.sh})"

cd $CONDA_PREFIX

RELEASE_URL="https://github.com/Rosemeis/pcangsd/archive/refs/tags/v.0.99.tar.gz"
wget -O- $RELEASE_URL | tar -xvzf-
cd pcangsd-v.0.99

python setup.py build_ext --inplace
#pip3 install -e .
pip3 install .

sed -i '1s/^/#!\/usr\/bin\/env python\n/' pcangsd.py
chmod +x pcangsd.py
cd $CONDA_PREFIX/bin && ln -s ../pcangsd-v.0.99/pcangsd.py ./pcangsd

#echo '#!/usr/bin/env python' > $CONDA_PREFIX/bin/pcangsd
#cat pcangsd.py >> $CONDA_PREFIX/bin/pcangsd
#chmod +x $CONDA_PREFIX/bin/pcangsd

#cd $CONDA_PREFIX && rm -rf pcangsd-v.0.99
