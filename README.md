# acoustics-kand
Kandidatarbete i beräkningsvetenskap och akustik


INSTRUKTIONER GIT:
*_______________________________________________________*
| GÖR ALLTID GIT PULL INNAN DU PUSHAR/ADDAR/COMMITTAR!! |
*-------------------------------------------------------*

--- Hämta uppdateringar --- 
1. git pull

--- Lägg till ändringar ---
1. git add .     ( . = allt, du kan också välja enskilda filer)
2. git commit -m "meddelande, t.ex. vad du ändrat"
3. git push

--- Kolla status ---
(säger om du har ändringar du inte har addat/committat)
1. git status

--- Skriv över dina ej committade ändringar ---
(dina ändringar försvinner inte helt, sparas bara på ett annat ställe)
1. git pull
2. (Du får ett meddelande om att du måste göra stash eller commit för dina ändringar)
3. git stash 

--- Skriv över din lokala mapp med main --- 
(om du råkat committa dina ändringar och det står att du är ahead of main branch)
1. git reset --hard origin/main



INSTRUKTIONER FÖR UPPMAX

--- Logga in på Uppmax ---
1. Öppna 2 terminaler
2. Gå till acoustics-kand i den ena (lokala terminalen)
3. I den andra: ssh alvafah@rackham.uppmax.uu.se
4. Skriv in ditt lösenord
5. Gå till vår mapp: cd /proj/snic2022-22-185/private/wallsource

--- Föra över filer att köra ---
1. I lokala terminalen: 
   scp FILNAMN alvafah@rackham.uppmax.uu.se:/proj/snic2022-22-185/private/wallsource
2. Skriv lösenord
3. Kolla i uppmax-terminalen om den kommit fram: 
   ls

--- Skicka jobb till Uppmax ---
1. För över dina filer till Uppmax (ovan)
2. Redigera batchfilen/skapa en ny batchfil som kör det/de program du vill köra
    2a. vim batch.sh
    2b. i           (insert mode)
    2c. ~ gör dina ändringar ~
    2d. esc         (backar ut från insert mode)
    2e. :w          (sparar)
    2f. :q          (stänger)
3. Kör batchfilen (här heter den batch.sh):
   sbatch batch.sh
   Du får ett jobbnummer, spara det
4. Kolla status på ditt jobb:
   jobinfo -u alvafah
5. När jobbet är klart skapas en fil som heter
   slurm-JOBBNUMMER.out
6. Öppna den för att se programmets utskrifter:
   cat slurm-JOBBNUMMER.out

--- Skicka tillbaka filer till datorn ---
1. Detta görs i lokala terminalen!!!
2. För att kopiera filen till där du är i lokala terminalen:
   scp alvafah@rackham.uppmax.uu.se:/proj/snic2022-22-185/private/wallsource/FILNAMN .
3. Vill du lägga den någon annanstans, byt ut punkten mot rätt adress

