FEMAG TS Postcalculation
************************

\femagCommand{import femagtools\_ts as ftts}

FEMAG-TS erstellt für die vtu-Files ein Verzeichnis mit dem Modellname und angegängtem "\_results\_x".
Beim Initialisieren des ts\_vtu-Moduls muss der Modellname und das Verzeichnis der vtu-Files angegeben werden.
\vspace{3mm} \\
\femagParameter{modelname = "model"} \\
\femagParameter{directory = modelname+"\_results\_1/"}
\smalllineskip
\femagCommand{vtu\_data = ftts.ts\_vtu(modelname,directory)}
