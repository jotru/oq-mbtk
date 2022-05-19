cd $ENV:USERPROFILE
$curDir = Get-Location
$workDir = $curDir.Path + "\mbtk"
Write-Host "Current Working Directory: $curDir"
Write-Host "Download and Installation of python embeddable package"
Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.8.10/python-3.8.10-embed-amd64.zip -OutFile .\python.zip
Expand-Archive -Path .\python.zip -DestinationPath .\mbtk -Verbose
Remove-Item .\python.zip
cd $workDir
Write-Host "Current Working Directory: $workDir"
Write-Host "Download and Installation of pip, wheel and setuptools"
Invoke-WebRequest -Uri https://bootstrap.pypa.io/get-pip.py -OutFile .\get-pip.py
$Env:PY_PIP = "$workDir\Scripts"
$Env:PY_LIBS = "$workDir\Lib;$workDir\Lib\site-package"
$Env:PY_HOME = "$workDir"
$Env:PATH="$workDir\Scripts;$ENV:PATH"
$Env:PYTHONUTF8=1
Set-Alias -Name python -Value $Env:PY_HOME\python.exe
Set-Alias -Name pip -Value $Env:PY_PIP\pip.exe
#
python .\get-pip.py
Move-Item .\python38._pth .\python38._pth.old
#
pip install pytest
Write-Host "clone of the branch $branch for oq-engine and install in developer mode"
git clone --depth=1 https://github.com/gem/oq-engine.git
cd .\oq-engine\
pip install -r .\requirements-py38-win64.txt
pip install -e .
cd ..
Write-Host "clone of the branch $branch for oq-mbtk and install in developer mode"
git clone --depth=1 https://github.com/GEMScienceTools/oq-mbtk.git
cd .\oq-mbtk\
pip install -r .\requirements_win64.txt
pip install -e .
Write-Host "End of installation"
Write-Host "Creation of symlink for bat files on the Desktop of user $ENV:USERNAME"
cd windows
Copy-Item  -Path .\oq-console.cmd -Destination "$ENV:USERPROFILE\Desktop"
Copy-Item  -Path .\oq-server.cmd -Destination "$ENV:USERPROFILE\Desktop"
cd $workDir
