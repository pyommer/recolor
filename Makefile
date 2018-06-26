CC = g++
CFLAGS = -g -Wall
LFLAGS = -g


ifeq ("$(shell uname)", "Darwin")
    LDFLAGS     = -framework Foundation -framework GLUT -framework OpenGL -lOpenImageIO -lm
else
    ifeq ("$(shell uname)", "Linux")
        LDFLAGS     = -L /usr/lib64/ -lglut -lGL -lGLU -lOpenImageIO -lm
    endif
endif

PROJECT = recolor
OBJS = recolor.o

#this generically compiles each .cpp to a .o file
%.o: %.cpp
    ${CC} -c ${CFLAGS} $<

#this does the linking step
all: recolor

recolor: recolor.o
    ${CC} ${LFLAGS} -o recolor recolor.o ${LDFLAGS}

#this will clean up all temporary files created by make all
clean:
    rm -f core.* *.o *~ ${PROJECT}
