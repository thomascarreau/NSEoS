CC = gcc
RM = rm
OBJ = etf_num.c
APP = etf_num

$(APP) : $(OBJ)
	$(CC) -o $(APP) $(OBJ) -lgsl -lgslcblas -lm

clean:
	$(RM) $(APP)
