

def loadSeq(fileName):
    name = None
    contents = []
    
    with open (fileName, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                name = line
            else:
                contents.append(line)
        
        return ("".join(contents))
       

def main():
    print(len(loadSeq("X73525.fa")))

    

if __name__ == "__main__":
    main()