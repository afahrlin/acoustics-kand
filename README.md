# acoustics-kand
Kandidatarbete i beräkningsvetenskap och akustik


Instruktioner för git-användning:

GÖR ALLTID GIT PULL INNAN DU PUSHAR/ADDAR/COMMITTAR!!

För att hämta uppdateringar:

1. git pull

För att lägga till ändringar:

1. git add .     ( . = allt, du kan också välja enskilda filer)

2. git commit -m "meddelande, t.ex. vad du ändrat"

3. git push

För att kolla status (säger om du har ändringar du inte har addat/committat):

1. git status

Om du vill kopiera  det som ligger i gitmolnet, så att dina egna lokala ändringar "försvinner"
(de försvinner inte helt, sparas bara på ett annat ställe)

1. git pull

2. (Du får ett meddelande om att du måste göra stash eller commit för dina ändringar)

3. git stash 

Ett annat sätt, om du råkat committa dina ändringar och det står att du är ahead of main branch:

1. git reset --hard origin/main
