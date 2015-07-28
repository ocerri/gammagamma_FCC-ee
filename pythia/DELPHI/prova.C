#include <dos.h>
#include <stdio.h>
#include <conio.h>

int main()
{
   printf("Message 1\n");
   sleep(2); //Parameter in sleep is in seconds
   printf("Message 2 a two seconds after Message 1");
   getch();
   return 0;
}
