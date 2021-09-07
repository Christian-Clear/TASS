# TAME

Term Analysis Made Easy (TAME) is designed to make the process of term analysis more user-friendly and promote the standardised storage and transfer of project files.
        
TAME uses a new implementation of the STRANS code to find matching lines in a linelist from a list of user-supplied energy levels. These matched lines are passed to the least-squares optimisation program (LOPT) to calculate optimised energy levels.
        
To download the latest TAME version or submit bug reports and improvement suggestions, please visit the project homepage. 

## Installation

The program is unstable and experimental right now, so installation is a
little tricky. Here's what I'm doing to install.

```bash
pip install -U .
sudo ln -s share/applications/com.codemouse92.thataway.desktop /usr/share/applications/
xdg-settings set default-web-browser thataway.desktop
```

(Please don't `sudo pip install`. Ever. Even here.)


