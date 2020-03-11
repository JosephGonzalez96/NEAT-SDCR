## NEAT implementation
This code is an implementation of the NEAT algorithm proposed by Stanley [in his paper](http://nn.cs.utexas.edu/downloads/papers/stanley.ec02.pdf). This algorithm uses the NEAT algorithm to control an Arduino robot to follow green colour objects. The original algorithm uses sexual reproduction  (**neat_sexual.c** file); this work proposes asexual reproduction (**neat_asexual.c** file) as a better alternative. The document that explains asexual reproduction and other innovations is available on [this website](https://sites.google.com/site/degreethesislorenaguachi/2020-joseph-gonzalez-self-driving-mini-robot-using-neat-algorithm). 


You can compile the **neat_sexual.c** file using the command 
```bash
gcc neat_sexual.c -o neat -lm
```
and the **neat_asexual.c** file using the command
```bash
gcc neat_asexual.c -o neat -lm
```
After training, the genome.txt and the resume.txt files will be generated; these files show a summary of the best neural networks in each generation and their structures. The asexual model usually reaches an accuracy of 95% in 5 minutes; nevertheless, it is important to adjust the hyper-parameters *population* and *max_generations* to reach an accuracy of 100% if it is possible. After training, the neural network has to be adapted to solve problems in big images; this process is not always simple. The neural network has to be adapted and place in the **side_computer.cpp** file, after that, it is necessary to install *eigen*, *bluez* and *opencv*, compile the **side_computer.cpp** program using
```bash
g++ -I /path/to/eigen/ -o computer side_computer.cpp -lbluetooth `pkg-config --cflags --libs opencv`
```
**side_computer.cpp** program will send instructions of movements by Bluetooth, and the Arduino program  **side_robot.ino** will receive these instructions to control the movement of the mini-robot