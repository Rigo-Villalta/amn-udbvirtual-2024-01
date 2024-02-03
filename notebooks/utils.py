from rich.console import Console
from rich.table import Table


def imprimir_tabla(lista):
    """
    Función que toma una lista de listas y retorna una tabla rich, 
    la primer lista tiene los encabezados de las columnas y las siguientes
    los datos.
    """

    # Reinicio la tabla y la consola, ya que si estan en el mismo namespace
    # va acumulando las tablas
    tabla = Table()
    for column in lista[0]:
        tabla.add_column(column)
    
    for row in lista[1:]:
        tabla.add_row(*row)
    
    console = Console()
    console.print(tabla)
    
