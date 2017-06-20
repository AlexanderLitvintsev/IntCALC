object Form6: TForm6
  Left = 367
  Top = 141
  BorderIcons = [biSystemMenu]
  BorderStyle = bsSingle
  Caption = 'База Элементов'
  ClientHeight = 161
  ClientWidth = 509
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -10
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnClose = FormClose
  OnCreate = FormCreate
  OnShow = FormShow
  PixelsPerInch = 96
  TextHeight = 13
  object DBGrid1: TDBGrid
    Left = 0
    Top = 0
    Width = 505
    Height = 121
    DataSource = DataSource1
    TabOrder = 0
    TitleFont.Charset = DEFAULT_CHARSET
    TitleFont.Color = clWindowText
    TitleFont.Height = -10
    TitleFont.Name = 'MS Sans Serif'
    TitleFont.Style = []
    Columns = <
      item
        Expanded = False
        FieldName = 'ElemName'
        Title.Caption = 'Наименование'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P1'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P2'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P3'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P4'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P5'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P6'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P7'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P8'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P9'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P10'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P11'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P12'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P13'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P14'
        Visible = True
      end
      item
        Expanded = False
        FieldName = 'P15'
        Visible = True
      end>
  end
  object Button1: TButton
    Left = 208
    Top = 128
    Width = 75
    Height = 25
    Caption = 'Выбрать'
    TabOrder = 1
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 296
    Top = 128
    Width = 75
    Height = 25
    Caption = 'Отмена'
    TabOrder = 2
    OnClick = Button2Click
  end
  object DBNavigator1: TDBNavigator
    Left = 0
    Top = 128
    Width = 192
    Height = 25
    DataSource = DataSource1
    VisibleButtons = [nbFirst, nbPrior, nbNext, nbLast, nbInsert, nbDelete]
    TabOrder = 3
  end
  object Table1: TTable
    DatabaseName = 'ElemBase'
    IndexFieldNames = 'No'
    TableName = 'elembase.db'
    Left = 16
    Top = 120
  end
  object DataSource1: TDataSource
    DataSet = Table1
    Left = 40
    Top = 120
  end
end
