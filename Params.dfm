object Form2: TForm2
  Left = 290
  Top = 233
  Anchors = []
  BorderIcons = []
  BorderStyle = bsSingle
  Caption = 'Параметры Элемента'
  ClientHeight = 224
  ClientWidth = 530
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -10
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  Position = poDesktopCenter
  OnActivate = FormActivate
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 392
    Top = 35
    Width = 32
    Height = 13
    Anchors = []
    Caption = 'Label1'
    Visible = False
  end
  object Panel2: TPanel
    Left = 178
    Top = 0
    Width = 178
    Height = 226
    TabOrder = 0
    object ParamsCheck: TCheckBox
      Left = 8
      Top = 8
      Width = 161
      Height = 17
      Caption = 'Приведеные параметры'
      TabOrder = 0
      OnClick = ParamsCheckClick
    end
  end
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 177
    Height = 225
    TabOrder = 1
  end
  object Button1: TButton
    Left = 372
    Top = 4
    Width = 61
    Height = 20
    Caption = 'ОК'
    TabOrder = 2
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 372
    Top = 28
    Width = 61
    Height = 20
    Caption = 'Отмена'
    TabOrder = 3
    OnClick = Button2Click
  end
end
