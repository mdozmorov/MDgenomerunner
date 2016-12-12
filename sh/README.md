# Misc scripts and notes for Mac/Linux shell

- `sqzspaces.sh` - removes spaces from file names. Use as `for file in *.pdf; do ./sqzspaces.sh "$file"; done`



# Misc notes

## Convert PDF to PNG

From [http://www.ademcan.net/?d=2013/04/10/14/23/32-how-to-convert-pdf-to-png-from-the-command-line-on-a-mac](http://www.ademcan.net/?d=2013/04/10/14/23/32-how-to-convert-pdf-to-png-from-the-command-line-on-a-mac)

`for file in *.pdf; do sips -s format png "$file" --out `basename $file .pdf`.png; done`

