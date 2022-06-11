import sys 
import os

def replaceInFile(fileName, destFileName, templateName : str, cloneName: str):
   with open(fileName, "r+") as text_file:
        texts = text_file.read()
        texts = texts.replace(templateName, cloneName)
        texts = texts.replace(templateName.capitalize(), cloneName.capitalize())
        texts = texts.replace(templateName.upper(), cloneName.upper())
   with open(destFileName, "w") as text_file:
        text_file.write(texts)


def copyAllFiles(solversDir, templateName, cloneName):
    # get list of files in dir
    baseDir = os.path.join(solversDir, templateName)
    baseFiles = [f for f in os.listdir(baseDir) if os.path.isfile(os.path.join(baseDir, f))]
    
    # create new dir
    cloneDir = os.path.join(solversDir, cloneName)
    if os.path.isdir(cloneDir):
      print(f"ERROR: {cloneDir} is an existing directory")
      exit(1)
    os.mkdir(cloneDir)

    # copy all files renaming each occurence of BASENAME
    for f in baseFiles:
      newFname = os.path.join(cloneDir, f.replace(templateName, cloneName))
      replaceInFile(os.path.join(baseDir, f), newFname, templateName, cloneName)
      print(f"Created {newFname}")
    return cloneDir


# execute function
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Syntax: createDriver.py inputtemplate outputname")
        exit(1)

    solversDir = os.path.abspath(os.path.dirname(__file__))
    baseName = sys.argv[1]
    newName = sys.argv[2]
    print(f"Cloning {os.path.join(solversDir, baseName)} into {os.path.join(solversDir, newName)}")

    destDir = copyAllFiles(solversDir, baseName, newName)

    print(f"\n\nNew driver files created in {destDir}")
    print(f"\nAdd the following line to ${os.path.join(solversDir, 'CMakeLists.txt')}:\n  add_ampl_backend({newName})")
    print("Solver libraries to link with are probably needed, find the relevant section in the documentation")