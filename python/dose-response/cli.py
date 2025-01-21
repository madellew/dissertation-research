import typer
from dissertation import dose_response_optim_boot
from typing_extensions import Annotated
app = typer.Typer()


@app.command()
def dose_response(
        models_to_fit: Annotated[
            list[dose_response_optim_boot.Model],
            typer.Argument(help="A model to fit. Multiple can be list at once")
        ]
):
    dose_response_optim_boot.main(models_to_fit=models_to_fit)

if __name__ == '__main__':
    app()
